/*
 * @file unphysical_masses.C
 * @brief create different plots to check what changed before and after VertexFit
 * @author Jennifer Puetz (j.puetz@fz-juelich.de)
 * @date 2015
 */

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "../common_jenny.cpp"

void unphysical_masses_boxgen(bool saveoutput=true, bool close=false, bool daughter=false){

	//****Input
	TString inPath = "/home/ikp1/puetz/panda/mysimulations/test/boxgenerator/lambda0/10000_events/";
	TString outPath  = inPath + "plots/";

	TFile * inFile = new TFile(inPath + "output_ana.root", "READ");
	TTree * ntpLambda0 = (TTree*) inFile->Get("ntpLambda0");

	TH1D * h_m_moth = new TH1D ("h_m_moth", "Mass of Lambda0 after Fit; m/ GeV/c^{2}; counts", 100, 1,1.22);
	ntpLambda0->Project("h_m_moth", "VtxFit_m", "McTruthMatch==1 && VtxFit_prob>0.01");

	TH1D * h_mreco_moth = new TH1D ("h_mreco_moth", "Mass of Lambda0 before Fit; m/ GeV/c^{2}; counts", 100, 1.06,1.3);
	ntpLambda0->Project("h_mreco_moth", "Lambda0_m", "McTruthMatch==1 && VtxFit_prob>0.01");

	TH1D * h_mdiff_moth = new TH1D ("h_mdiff_moth", "Mass diff of Lambda0; m_{reco}-m_{fit}/ GeV/c^{2}; counts", 100, -0.2,0.2);
	ntpLambda0->Project("h_mdiff_moth", "Lambda0_m-VtxFit_m", "McTruthMatch==1 && VtxFit_prob>0.1");

	TH1D * h_p_moth = new TH1D ("h_p_moth", "Momentum of Lambda0 after Fit; m/ GeV/c^{2}; counts", 100, 0,5);
	ntpLambda0->Project("h_p_moth", "VtxFit_p", "McTruthMatch==1 && VtxFit_prob>0.1");

	TH1D * h_preco_moth = new TH1D ("h_preco_moth", "Mometum of Lambda0 before Fit; m/ GeV/c^{2}; counts", 100, 0,5);
	ntpLambda0->Project("h_preco_moth", "Lambda0_p", "McTruthMatch==1 && VtxFit_prob>0.1");

	TH1D * h_pdiff_moth = new TH1D ("h_pdiff_moth", "Momentum diff of Lambda0; p_{reco}-p_{fit}/ GeV/c; counts", 100, -0.2,0.2);
	ntpLambda0->Project("h_pdiff_moth", "Lambda0_p-VtxFit_p", "McTruthMatch==1 && VtxFit_prob>0.1");

	TH1D * h_e_moth = new TH1D ("h_e_moth", "Energy of Lambda0 after Fit; m/ GeV/c^{2}; counts", 100, 0,5);
	ntpLambda0->Project("h_e_moth", "VtxFit_e", "McTruthMatch==1 && VtxFit_prob>0.1");

	TH1D * h_ereco_moth = new TH1D ("h_ereco_moth", "Energy of Lambda0 before Fit; m/ GeV/c^{2}; counts", 100, 0,5);
	ntpLambda0->Project("h_ereco_moth", "Lambda0_e", "McTruthMatch==1 && VtxFit_prob>0.1");

	TH1D * h_ediff_moth = new TH1D ("h_ediff_moth", "Energy diff of Lambda0; e_{reco}-e_{fit}/ GeV; counts", 100, -0.2,0.2);
	ntpLambda0->Project("h_ediff_moth", "Lambda0_e-VtxFit_e", "McTruthMatch==1 && VtxFit_prob>0.1");

//	jenny::CreateDrawAndSaveHistogram(h_m_moth, outPath, "h_mreco_moth", saveoutput, close);
	jenny::CreateDrawAndSaveNHistograms(h_m_moth, h_mreco_moth, "after fit", "before fit", outPath, "h_m_before_and_after_fit", saveoutput, close);
	jenny::CreateDrawAndSaveHistogram(h_mdiff_moth, outPath, "h_mdiff_moth", saveoutput, close);
	jenny::CreateDrawAndSaveNHistograms(h_preco_moth, h_p_moth, "before fit", "after fit", outPath, "h_p_before_and_after_fit", saveoutput, close);
	jenny::CreateDrawAndSaveHistogram(h_pdiff_moth, outPath, "h_pdiff_moth", saveoutput, close);
	jenny::CreateDrawAndSaveNHistograms(h_ereco_moth, h_e_moth, "before fit", "after fit", outPath, "h_e_before_and_after_fit", saveoutput, close);
	jenny::CreateDrawAndSaveHistogram(h_ediff_moth, outPath, "h_ediff_moth", saveoutput, close);

	if (daughter){


		TH1D * h_mdiff_dau0 = new TH1D ("h_mdiff_dau0", "Mass diff of Daughter 1; m_{reco}-m_{fit}/ GeV/c^{2}; counts", 100, -0.05,0.05);
		ntpLambda0->Project("h_mdiff_dau0", "Lambda0_d0m-VtxFit_d0m", "McTruthMatch==1");

		TH1D * h_pdiff_dau0 = new TH1D ("h_pdiff_dau0", "Momentum diff of Daughter 1; p_{reco}-p_{fit}/ GeV/c; counts", 100, -0.05,0.05);
		ntpLambda0->Project("h_pdiff_dau0", "Lambda0_d0p-VtxFit_d0p", "McTruthMatch==1");

		TH1D * h_ediff_dau0 = new TH1D ("h_ediff_dau0", "Energy diff of Daughter 1; e_{reco}-e_{fit}/ GeV; counts", 100, -0.05,0.05);
		ntpLambda0->Project("h_ediff_dau0", "Lambda0_d0e-VtxFit_d0e", "McTruthMatch==1");

		jenny::CreateDrawAndSaveHistogram(h_mdiff_dau0, outPath, "h_mdiff_dau0", saveoutput, close);
		jenny::CreateDrawAndSaveHistogram(h_pdiff_dau0, outPath, "h_pdiff_dau0", saveoutput, close);
		jenny::CreateDrawAndSaveHistogram(h_ediff_dau0, outPath, "h_ediff_dau0", saveoutput, close);




		TH1D * h_mdiff_dau1 = new TH1D ("h_mdiff_dau1", "Mass diff of Daughter 2; m_{reco}-m_{fit}/ GeV/c^{2}; counts", 100, -0.05,0.05);
		ntpLambda0->Project("h_mdiff_dau1", "Lambda0_d1m-VtxFit_d1m", "McTruthMatch==1");

		TH1D * h_pdiff_dau1 = new TH1D ("h_pdiff_dau1", "Momentum diff of Daughter 2; p_{reco}-p_{fit}/ GeV/c; counts", 100, -0.05,0.05);
		ntpLambda0->Project("h_pdiff_dau1", "Lambda0_d1p-VtxFit_d1p", "McTruthMatch==1");

		TH1D * h_ediff_dau1 = new TH1D ("h_ediff_dau1", "Energy diff of Daughter 2; e_{reco}-e_{fit}/ GeV; counts", 100, -0.05,0.05);
		ntpLambda0->Project("h_ediff_dau1", "Lambda0_d1e-VtxFit_d1e", "McTruthMatch==1");

		jenny::CreateDrawAndSaveHistogram(h_mdiff_dau1, outPath, "h_mdiff_dau1", saveoutput, close);
		jenny::CreateDrawAndSaveHistogram(h_pdiff_dau1, outPath, "h_pdiff_dau1", saveoutput, close);
		jenny::CreateDrawAndSaveHistogram(h_ediff_dau1, outPath, "h_ediff_dau1", saveoutput, close);


	}

	if(close) exit(0);

}
