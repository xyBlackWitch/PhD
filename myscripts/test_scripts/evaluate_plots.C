/**@file evaluate_plots.C
 * @author Jennifer Puetz (j.puetz@fz-juelich.de)
 * @date 2016
 * @brief evaluate comparison plots
 * @details This macro evaluates two different histograms from root files for comparison
 *
 */

#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "../common_jenny.cpp"


void evaluate_plots(TString file1="", TString name1="", TString file2="", TString name2=""){

	//*** Create Canvas
	TCanvas * c = new TCanvas("c", "compare histograms", 0,0,800,500);

	//*** Create new Legend
	TLegend * legend = new TLegend(0.68,0.74,0.99,0.88, "");

	//*** Get files
	TFile * f1 = new TFile(file1, "READ");
	TFile * f2 = new TFile(file2, "READ");

	//*** Get histograms
	TH1D * h1 = f1->Get("h_eff");
	TH1D * h2 = f2->Get("h_eff");

	//*** Draw histograms
	jenny::comparisonStyle(h1);
	h1->SetLineColor(kBlue);
	h1->SetMarkerColor(kBlue);
	legend->AddEntry(h1, name1, "pl");

	h1->Draw("PL");

	h2->SetLineColor(kRed);
	h2->SetMarkerColor(kRed);
	legend->AddEntry(h2, name2, "pl");

	h2->Draw("SAME PL");

	legend->Draw();

	c->Print("comparison_PID.root");
	c->Print("comparison_PID.pdf");


}
