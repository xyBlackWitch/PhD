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
#include "common_jenny.cpp"

void vertex_studies_XiMinus(int nevts=0, bool saveoutput=true, bool close=false){

	//Paths

	TString inPath = "/private/puetz/mysimulations/test/boxgenerator/Xi/run2/old/";
	TString outPath = inPath + "plots/";


	//Input File

	TString inputFile = inPath + "output_ana_oldidealtracking_new.root";

	//open input
	TFile * file = new TFile(inputFile, "READ");

	TTree * ntpXiMinus = (TTree*) file->Get("ntpXiMinus");
	TTree * ntpMC = (TTree*) file->Get("ntpMC");

	//****Projection for XiMinus

	// simulated vertex position for XiMinus

	TH1D* h_vx_mc = new TH1D("h_vx_mc", "x coordinate of the vertex position; x[cm]; counts", 50,-5,5);
	ntpMC->Project("h_vx_mc", "x", "moth==0 && pdg==3312");

	TH1D* h_vy_mc = new TH1D("h_vy_mc", "y coordinate of the vertex position; y[cm]; counts", 50,-5,5);
	ntpMC->Project("h_vy_mc", "y", "moth==0 && pdg==3312");

	TH1D* h_vz_mc = new TH1D("h_vz_mc", "z coordinate of the vertex position; z[cm]; counts", 50,-5,5);
	ntpMC->Project("h_vz_mc", "z", "moth==0 && pdg==3312");



	//reconstructed vertex position for XiMinus
	TH1D* h_vx_reco = new TH1D("h_vx_reco", "x coordinate of the vertex position; x[cm]; counts", 50,-5,5);
	ntpXiMinus->Project("h_vx_reco", "XiMinusFit_vx", "McTruthMatch==1");

	TH1D* h_vy_reco = new TH1D("h_vy_reco", "y coordinate of the vertex position; y[cm]; counts", 50,-5,5);
	ntpXiMinus->Project("h_vy_reco", "XiMinusFit_vy", "McTruthMatch==1");

	TH1D* h_vz_reco = new TH1D("h_vz_reco", "z coordinate of the vertex position; z[cm]; counts", 50,-5,5);
	ntpXiMinus->Project("h_vz_reco", "XiMinusFit_vz", "McTruthMatch==1 ");



	//resolution of reconstructed vertex position

	TH1D* h_vx_res = new TH1D("h_vx_res", "x resolution of the vertex position; #Delta_{x}[cm]; counts", 100,-1,1);
	ntpXiMinus->Project("h_vx_res", "XiMinusFit_vx-MCTruth_vx", "McTruthMatch==1 ");

	TH1D* h_vy_res = new TH1D("h_vy_res", "y resolution of the vertex position; #Delta_{y}[cm]; counts", 100,-1,1);
	ntpXiMinus->Project("h_vy_res", "XiMinusFit_vy-MCTruth_vy", "McTruthMatch==1 ");

	TH1D* h_vz_res = new TH1D("h_vz_res", "z resolution of the vertex position; #Delta_{z}[cm]; counts", 100,-1.50,1.50);
	ntpXiMinus->Project("h_vz_res", "XiMinusFit_vz-MCTruth_vz", "McTruthMatch==1");


	//resolution of reconstructed vertex position and the poca

	TH1D* h_vx_poc_res = new TH1D("h_vx_poc_res", "differenz reconstructed vertex position and poca position; x_{reco} - x_{poca}[cm]; counts", 100,-1,1);
	ntpXiMinus->Project("h_vx_poc_res", "XiMinusFit_vx - XiMinus_pocvx", "McTruthMatch==1 ");

	TH1D* h_vy_poc_res = new TH1D("h_vy_poc_res", "differenz reconstructed vertex position and poca position; y_{reco} - y_{poca}[cm]; counts", 100,-1,1);
	ntpXiMinus->Project("h_vy_poc_res", "XiMinusFit_vy - XiMinus_pocvy", "McTruthMatch==1 ");

	TH1D* h_vz_poc_res = new TH1D("h_vz_poc_res", "differenz reconstructed vertex position and poca position; z_{reco} - z_{poca}[cm]; counts", 100,-1.5,1.5);
	ntpXiMinus->Project("h_vz_poc_res", "XiMinusFit_vz - XiMinus_pocvz", "McTruthMatch==1");




	//ratio plot of vertex postion for XiMinus
	TH1D* h_reco_div_mc_x = new TH1D("h_reco_div_mc_x", "ratio reco/MC for x position of vertex; x[cm]; ratio reco/MC", 50,-5,5);
	h_reco_div_mc_x->Divide(h_vx_reco, h_vx_mc);
	h_reco_div_mc_x->SetStats(0);

	TH1D* h_reco_div_mc_y = new TH1D("h_reco_div_mc_y", "ratio reco/MC for y position of vertex; y[cm]; ratio reco/MC", 50,-5,5);
	h_reco_div_mc_y->Divide(h_vy_reco, h_vy_mc);
	h_reco_div_mc_y->SetStats(0);

	TH1D* h_reco_div_mc_z = new TH1D("h_reco_div_mc_z", "ratio reco/MC for z position of vertex; z[cm]; ratio reco/MC", 50,-5,5);
	h_reco_div_mc_z->Divide(h_vz_reco, h_vz_mc);
	h_reco_div_mc_z->SetStats(0);



	//2-dimensional plots vertex

	TH2D * h_vxy_vz_reco = new TH2D("h_vxy_vz_reco", "Vertex ; z[cm]; xy[cm]", 50,0,5,50,0,5);
	ntpXiMinus->Project("h_vxy_vz_reco", "sqrt(XiMinusFit_vx**2+XiMinusFit_vy**2):XiMinusFit_vz", "McTruthMatch==1 ");



	//Goodness of fit

	TH1D * h_chi2 = new TH1D("h_chi2", "#chi^{2} distribution for vertex Fit; #chi^{2}; counts", 50,0,200);
	ntpXiMinus->Project("h_chi2", "FitVertex_Chi2", "McTruthMatch==1");

	TH1D * h_Prob = new TH1D("h_Prob", "Probability for vertex Fit; Prob; counts", 50,0,1);
	ntpXiMinus->Project("h_Prob", "FitVertex_Prob", "McTruthMatch==1");



	//****Create Canvas, draw histogram and save it

	jenny::CreateDrawAndSaveHistogram(h_vx_mc, outPath, "h_vx_mc", saveoutput, close);
	jenny::CreateDrawAndSaveHistogram(h_vy_mc, outPath, "h_vy_mc", saveoutput, close);
	jenny::CreateDrawAndSaveHistogram(h_vz_mc, outPath, "h_vz_mc", saveoutput, close);

	jenny::CreateDrawAndSaveHistogram(h_vx_reco, outPath, "h_vx_reco", saveoutput, close);
	jenny::CreateDrawAndSaveHistogram(h_vy_reco, outPath, "h_vy_reco", saveoutput, close);
	jenny::CreateDrawAndSaveHistogram(h_vz_reco, outPath, "h_vz_reco", saveoutput, close);

	bool noAutorange =false;

	jenny::CreateDrawAndSaveHistogramWithFit(h_vx_res, outPath, "h_vx_res_fit", saveoutput, close, noAutorange, 0.09, 1);
	jenny::CreateDrawAndSaveHistogramWithFit(h_vy_res, outPath, "h_vy_res_fit", saveoutput, close, noAutorange, 0.08, 1);
	jenny::CreateDrawAndSaveHistogramWithFit(h_vz_res, outPath, "h_vz_res_fit", saveoutput, close, noAutorange, 0.15,1);

//	jenny::CreateDrawAndSaveHistogram(h_vx_res, outPath, "h_vx_res", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_vy_res, outPath, "h_vy_res", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_vz_res, outPath, "h_vz_res", saveoutput, close);

	jenny::CreateDrawAndSaveHistogramWithFit(h_vx_poc_res, outPath, "h_vx_poc_res", saveoutput, close);//, noAutorange, 0.16, 1);
	jenny::CreateDrawAndSaveHistogramWithFit(h_vy_poc_res, outPath, "h_vy_poc_res", saveoutput, close, noAutorange, 0.02, 2);
	jenny::CreateDrawAndSaveHistogramWithFit(h_vz_poc_res, outPath, "h_vz_poc_res", saveoutput, close, noAutorange, 0.15,5);

//	jenny::CreateDrawAndSaveHistogram(h_reco_div_mc_x, outPath, "h_reco_div_mc_x", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_reco_div_mc_y, outPath, "h_reco_div_mc_y", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_reco_div_mc_z, outPath, "h_reco_div_mc_z", saveoutput, close);


	jenny::CreateDrawAndSaveHistogram(h_vxy_vz_reco, outPath, "h_vxy_vz_reco", saveoutput, close);

	jenny::CreateDrawAndSaveHistogram(h_chi2, outPath, "h_chi2", saveoutput, close);
	jenny::CreateDrawAndSaveHistogram(h_Prob, outPath, "h_Prob", saveoutput, close);

	if(close) exit(0);

}
