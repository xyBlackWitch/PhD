/***
*@file cut_optimization.C
*@author Jenny Puetz (j.puetz@fz-juelich.de)
*@date 2016
*@brief macro to find optimized cut
*@details This macro evaluates to optimum for the probability cut
*/

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1.h"

void cut_optimization(TString sig = "", TString bg = "", TString ntuple = "", int steps=100, double br=0.639){
	
	//*** input file
	TFile * sigData = new TFile(sig, "READ");
	TTree * signal = (TTree*) sigData->Get(ntuple);
	
	TFile * bkgData = new TFile(bg, "READ");
	TTree * bkg = (TTree*) bkgData->Get(ntuple);
	
	int nevts_sig = signal->GetEntriesFast();
	int nevts_bkg = bkg->GetEntriesFast();
	

	float factor = (float) nevts_sig*60e-3/(nevts_bkg*1e-6)*br;
	
	cout << "The scaling factor for background is " << factor << endl;
	
	gStyle->SetOptStat(0);

	TH1D * h_sig = new TH1D("h_sig", "probability; probability; counts", 100, 0, 1);
	TH1D * h_bkg = new TH1D("h_bkg", "probability; probability; counts", 100, 0, 1);
	TH1D * h_signif = new TH1D("h_signif", "probability; probability; N_{sig}/sqrt(N_{sig}+N_{bg}*B)", (steps-1), 0, 1);

	for (int i=1; i<steps; i++){

		double step = (double) i/100;
		TString cut = TString::Format("VtxFit_prob>%.2f", step);

		signal->Project("h_sig","VtxFit_prob", cut+ "& HitTag==1");
		int sig_evts = h_sig->GetEntries();

		bkg->Project("h_bkg","VtxFit_prob", cut+ "& Lambda0_HitTag==1");
		int bkg_evts = h_bkg->GetEntries();

		float signif = sig_evts/sqrt(sig_evts + bkg_evts*factor); //(bkg_evts*factor); //


		h_signif->Fill(step, signif);


	}
	TCanvas * c = new TCanvas("c", "significance", 0,0,800,500);
	
	int binx = h_signif->GetMaximumBin();
	float max = (float) binx/steps;
	float ymax =(int) TMath::Ceil(h_signif->GetMaximum())*1.1; // h_signif->GetMaximum() *1.1; //

	cout << "Maximum: (" << max << ", " << h_signif->GetMaximum() << ")" << endl;

	h_signif->GetYaxis()->SetRangeUser(0,ymax);

	TLine * line = new TLine(max, 0, max, ymax);
	h_signif->Draw("C");
	line->Draw("SAME");


	TCanvas * c1 = new TCanvas("c1", "prob", 0,0,800,500);
	c1->SetLogy();

	TLegend * legend = new TLegend(0.74,0.8,0.97,0.87, "");


	bkg->Project("h_bkg","VtxFit_prob");
	h_bkg->Scale(factor);
	h_bkg->GetYaxis()->SetRangeUser(1e1,1e10);
	h_bkg->SetLineColor(kBlue);
	h_bkg->SetFillColor(kBlue);
	h_bkg->SetFillStyle(3004);
	legend->AddEntry(h_bkg, "Background", "l");
	h_bkg->Draw();


	signal->Project("h_sig","VtxFit_prob");
	h_sig->SetLineColor(kGreen);
	h_sig->SetFillColor(kGreen);
	h_sig->SetFillStyle(3004);
	legend->AddEntry(h_sig, "Signal", "l");
	h_sig->Draw("SAME");

	legend->Draw();




}
