/**
 * @file vertex_studies_lambda0.C
 * @author Jennifer Puetz (j.puetz@fz-juelich.de)
 * @date 2015
 * @brief Creates histograms for vertex studies for Lambda0 - AntiLambda production
 * @details Methode creates Histogramms for vertex studies for Lambda0 - AntiLambda production. Data is comming from the analysis which is stored in a file called output_ana.root
 * Basically ROOT.
 */


class RhoTuple;

#include "common_jenny.cpp"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"


void vertex_studies_pbarp_to_2Lamda0(int nevts=0, bool saveoutput=true, bool close=false){

	//Paths
	TString inPath = "/private/puetz/mysimulations/analysis/pbarp_lambda0_antilambda0/1000_events/idealtracking/";
	TString outPath = inPath + "plots/";


	//Input file
	TFile * file = new TFile(inPath + "output_ana_oldidealtracking.root", "READ");

	//Reading data from file
	TTree * ntpMC = (TTree*) file->Get("ntpMC");
	TTree * ntpBeam = (TTree*) file->Get("ntpCrossCheck");

	//***Projection for the Beam

	//simulated vertex position for Beam

	TH1D * h_vx_mc = new TH1D("h_vx_mc", "x coordinate of vertex position; x[cm]; counts", 85,-5,5);
	ntpMC->Project("h_vx_mc", "x", "moth==0 && (pdg==3122 || pdg==-3122) ");

	TH1D * h_vy_mc = new TH1D("h_vy_mc", "y coordinate of vertex position; y[cm]; counts", 85,-5,5);
	ntpMC->Project("h_vy_mc", "y", "moth==0 && (pdg==3122 || pdg==-3122)");

	TH1D * h_vz_mc = new TH1D("h_vz_mc", "z coordinate of vertex position; z[cm]; counts", 85,0,5);
	ntpMC->Project("h_vz_mc", "z", "moth==0 && (pdg==3122 || pdg==-3122)");


	//reconstructed vertex position for Beam
	TH1D * h_vx_reco = new TH1D("h_vx_reco", "x coordinate of vertex position; x[cm]; counts", 58,-5,5);
	ntpBeam->Project("h_vx_reco", "VtxFit_vx", "McTruthMatch==1");

	TH1D * h_vy_reco = new TH1D("h_vy_reco", "y coordinate of vertex position; y[cm]; counts", 58,-5,5);
	ntpBeam->Project("h_vy_reco", "VtxFit_vy", "McTruthMatch==1");

	TH1D * h_vz_reco = new TH1D("h_vz_reco", "z coordinate of vertex position; z[cm]; counts", 58,-5,5);
	ntpBeam->Project("h_vz_reco", "VtxFit_vz", "McTruthMatch==1");


	//resolution vertex position for Beam
	TH1D * h_vx_res = new TH1D("h_vx_res", "x resolution of vertex position; #Delta_{x}[cm]; counts", 58,-5,5);
	ntpBeam->Project("h_vx_res", "VtxFit_mcdiffvx", "McTruthMatch==1");

	TH1D * h_vy_res = new TH1D("h_vy_res", "y resolution of vertex position; #Delta_{y}[cm]; counts", 58,-5,5);
	ntpBeam->Project("h_vy_res", "VtxFit_mcdiffvy", "McTruthMatch==1");

	TH1D * h_vz_res = new TH1D("h_vz_res", "z resolution of vertex position; #Delta_{z}[cm]; counts", 58,-5,5);
	ntpBeam->Project("h_vz_res", "VtxFit_mcdiffvz", "McTruthMatch==1");



	//2-dimensional plots for poca and vertex

	TH2D * h_vxy_vz_reco = new TH2D("h_vxy_vz_reco", "Vertex ; z[cm]; r[cm]", 58,0,5,58,0,5);
	ntpBeam->Project("h_vxy_vz_reco", "sqrt(VtxFit_vx**2+VtxFit_vy**2):VtxFit_vz", "McTruthMatch==1");


	//Goodness of fit

	TH1D * h_chi2 = new TH1D("h_chi2", "#chi^{2} distribution for vertex Fit; #chi^{2}; counts", 58,0,10);
	ntpBeam->Project("h_chi2", "VtxFit_chisq", "McTruthMatch==1");

	TH1D * h_Prob = new TH1D("h_Prob", "Probability for vertex Fit; Prob; counts", 58,0,1.1);
	ntpBeam->Project("h_Prob", "VtxFit_prob", "McTruthMatch==1");


	//*** create, draw and save histogram

	jenny::CreateDrawAndSaveHistogram(h_vx_mc, outPath, "h_vx_mc", saveoutput, close);
	jenny::CreateDrawAndSaveHistogram(h_vy_mc, outPath, "h_vy_mc", saveoutput, close);
	jenny::CreateDrawAndSaveHistogram(h_vz_mc, outPath, "h_vz_mc", saveoutput, close);

	jenny::CreateDrawAndSaveHistogram(h_vx_reco, outPath, "h_vx_reco", saveoutput, close);
	jenny::CreateDrawAndSaveHistogram(h_vy_reco, outPath, "h_vy_reco", saveoutput, close);
	jenny::CreateDrawAndSaveHistogram(h_vz_reco, outPath, "h_vz_reco", saveoutput, close);

	jenny::CreateDrawAndSaveHistogram(h_vx_res, outPath, "h_vx_res", saveoutput, close);
	jenny::CreateDrawAndSaveHistogram(h_vy_res, outPath, "h_vy_res", saveoutput, close);
	jenny::CreateDrawAndSaveHistogram(h_vz_res, outPath, "h_vz_res", saveoutput, close);


	jenny::CreateDrawAndSaveHistogram(h_vxy_vz_reco, outPath, "h_vxy_vz_reco", saveoutput, close);


	jenny::CreateDrawAndSaveHistogram(h_chi2, outPath, "h_chi2", saveoutput, close);
	jenny::CreateDrawAndSaveHistogram(h_Prob, outPath, "h_Prob", saveoutput, close);

	if(close) exit(0);

}
