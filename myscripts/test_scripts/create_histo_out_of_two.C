/**
 * @file create_histo_out_of_two.C
 * @author Jennifer Puetz (j.puetz@fz-juelich.de)
 * @date 2016
 * @brief Creates a histogram out of two histograms.
 * @details Creates a histogram with legend out of two histograms by reading two *.root-files.
 * Basically ROOT.
 */

void create_histo_out_of_two(TString infile1, TString histoname1, TString legend1,  TString infile2, TString histoname2, TString legend2){


	gStyle->SetOptStat(0);

	//************* Read files *****************
	TFile* inFile1 = new TFile(infile1, "READ");
	TFile* inFile2 = new TFile(infile2, "READ");


	//************* Create histograms *****************
	TH1D* histo1 = new TH1D("histo1", "reconstruction efficiencies; particle type; counts", 6, 0, 6);
//	TH1D* histo2 = new TH1D("histo2", "reconstruction efficiencies; particle type; counts", 6, 0, 6);
	TH1D* histo2 = (TH1D*) histo1->Clone("histo2");

	//************* Get entries from file *****************
	histo1 = (TH1D*) inFile1->Get(histoname1);
	histo1->SetLineColor(kBlue);

	histo1->GetXaxis()->SetLabelSize(0.06);
	histo1->GetXaxis()->SetTitleSize(0.05);

	histo1->GetYaxis()->SetTitleSize(0.05);

	histo2 = (TH1D*) inFile2->Get(histoname2);
	histo2->SetLineColor(kRed);

	//************* Create TCanvas and Legend **************
	TCanvas *c = new TCanvas("c", "reco eff", 0,0,1000,600);
	TLegend *l = new TLegend(0.751,0.748,0.978,0.869,"");

	l->AddEntry(histo1, legend1, "l");
	l->AddEntry(histo2, legend2, "l");

	histo1->Draw();
	histo2->Draw("SAME");
	l->Draw();
}
