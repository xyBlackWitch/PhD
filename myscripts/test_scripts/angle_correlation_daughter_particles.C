
void angle_correlation_daughter_particles(){

	TStopwatch timer;

	TString inPath ="";// "/home/ikp1/puetz/panda/mysimulations/test/boxgenerator/lambda0/10000_events/";

	TFile * inFile = new TFile(inPath + "output_ana.root", "READ");
	TTree * lambda0 = (TTree*) inFile->Get("fntpLambda0");

	gStyle->SetOptStat(0);

	TF1 * f = new TF1("f", "x", -1, 1);


	TH2D * h_cos_tht = new TH2D("h_cos_tht", "Correlation for tht of daughters; cos(tht_{d0});cos(tht_{d1})", 1000,-1,1, 1000,-1,1);
	lambda0->Project("h_cos_tht", "cos(Lambda0_d1tht):cos(Lambda0_d0tht)", "McTruthMatch==1");

	TH1D * h_tht = new TH1D("h_tht", "Correlation for tht of daughters; tht_{d0}-tht_{d1}/rad; counts", 100,-2,2);
	lambda0->Project("h_tht", "cos(Lambda0_d1tht)-cos(Lambda0_d0tht)", "McTruthMatch==1 && abs(Lambda0_d1phi-Lambda0_d0phi)<0.002");


	TH1D * h_phi = new TH1D("h_phi", "Correlation for phi of daughters; phi_{d0}-phi_{d1}/rad; counts", 100,-0.1,0.1);
	lambda0->Project("h_phi", "Lambda0_d1phi-Lambda0_d0phi", "McTruthMatch==1 && abs(Lambda0_d1phi-Lambda0_d0phi)<0.002");

	TCanvas * c_cos_tht = new TCanvas("c_cos_tht", "c_cos_tht", 0,0,800,500);
	h_cos_tht->Draw("COLZ");
	f->Draw("SAME");

	TCanvas * c_tht = new TCanvas("c_tht", "c_tht", 0,0,800,500);

	h_tht->Draw();



	TCanvas * c_phi = new TCanvas("c_phi", "c_phi", 0,0,800,500);

	h_phi->Draw();

}
