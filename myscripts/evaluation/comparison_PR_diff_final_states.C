#include <fstream>
#include <iostream>
#include "/home/ikp1/puetz/panda/PandaSoftware/pandaroot/trunk/source/macro/setPandaStyle.C"

void comparison_PR_diff_final_states(TString input_ideal = "", TString input_real=""){
	if(input_ideal=="" || input_real==""){
		cout << "Please insert two txt-file names" << endl;
	}

	std::ifstream filelist_ideal(input_ideal);

	if (!filelist_ideal.is_open()){
		cerr << "failed in open." << endl;
		exit();
	}

	TChain MC_ideal ("ntpMC");
	TChain PiMinus_ideal ("ntpPiMinus");
	TChain PiPlus_ideal ("ntpPiPlus");
	TChain KaonMinus_ideal ("ntpKaonMinus");
	TChain KaonPlus_ideal ("ntpKaonPlus");
	TChain Proton_ideal ("ntpProton");
	TChain AntiProton_ideal ("ntpAntiProton");


	char line_ideal[10000];

	while (!filelist_ideal.eof()){

		filelist_ideal.getline(line_ideal, sizeof(line_ideal));
		if (strcmp(line_ideal, "") == 0) continue;

		MC_ideal.Add(line_ideal);
		PiMinus_ideal.Add(line_ideal);
		PiPlus_ideal.Add(line_ideal);
		KaonMinus_ideal.Add(line_ideal);
		KaonPlus_ideal.Add(line_ideal);
		Proton_ideal.Add(line_ideal);
		AntiProton_ideal.Add(line_ideal);

	}

	std::ifstream filelist_real(input_real);

	if (!filelist_real.is_open()){
		cerr << "failed in open." << endl;
		exit();
	}

	TChain MC_real ("ntpMC");
	TChain PiMinus_real ("ntpPiMinus");
	TChain PiPlus_real ("ntpPiPlus");
	TChain KaonMinus_real ("ntpKaonMinus");
	TChain KaonPlus_real ("ntpKaonPlus");
	TChain Proton_real ("ntpProton");
	TChain AntiProton_real ("ntpAntiProton");


	char line_real[10000];

	while (!filelist_real.eof()){

		filelist_real.getline(line_real, sizeof(line_real));
		if (strcmp(line_real, "") == 0) continue;

		MC_real.Add(line_real);
		PiMinus_real.Add(line_real);
		PiPlus_real.Add(line_real);
		KaonMinus_real.Add(line_real);
		KaonPlus_real.Add(line_real);
		Proton_real.Add(line_real);
		AntiProton_real.Add(line_real);

	}

	setPandaStyle();


	double nevents_ideal = MC_ideal.GetEntries()/3;
	double nevents_real = MC_real.GetEntries()/3;

	cout << nevents_ideal << endl;
	cout << nevents_real << endl;

	TString cut = "Mother==88888 & McTruthMatch==1 & ";
	TString cut_real = "Mother==88888 & McTruthMatch==1";


	//Histograms
	// ideal
	TH1D * h_pim_ideal_tht = new TH1D("h_pim_ideal_tht", "tht", 200,0,3.5);
	PiMinus_ideal.Project("h_pim_ideal_tht", "piminus_tht", cut+"piminus_HitTag==1");
	double pim_ideal = h_pim_ideal_tht->GetEntries();

	TH1D * h_pip_ideal_tht = new TH1D("h_pip_ideal_tht", "tht", 200,0,3.5);
	PiPlus_ideal.Project("h_pip_ideal_tht", "piplus_tht", cut+"piplus_HitTag==1");
	double pip_ideal = h_pip_ideal_tht->GetEntries();

	TH1D * h_prot_ideal_tht = new TH1D("h_prot_ideal_tht", "tht", 200,0,3.5);
	Proton_ideal.Project("h_prot_ideal_tht", "proton_tht", cut+"proton_HitTag==1");
	double prot_ideal = h_prot_ideal_tht->GetEntries();

	TH1D * h_aprot_ideal_tht = new TH1D("h_aprot_ideal_tht", "tht", 200,0,3.5);
	AntiProton_ideal.Project("h_aprot_ideal_tht", "AntiProton_tht", cut+"AntiProton_HitTag==1");
	double aprot_ideal = h_aprot_ideal_tht->GetEntries();

	TH1D * h_km_ideal_tht = new TH1D("h_km_ideal_tht", "tht", 200,0,3.5);
	KaonMinus_ideal.Project("h_km_ideal_tht", "kaonminus_tht", cut+"kaonminus_HitTag==1");
	double km_ideal = h_km_ideal_tht->GetEntries();

	TH1D * h_kp_ideal_tht = new TH1D("h_kp_ideal_tht", "tht", 200,0,3.5);
	KaonPlus_ideal.Project("h_kp_ideal_tht", "kaonplus_tht", cut+"kaonplus_HitTag==1");
	double kp_ideal = h_kp_ideal_tht->GetEntries();


	//real
	TH1D * h_pim_real_tht = new TH1D("h_pim_real_tht", "tht", 200,0,3.5);
	PiMinus_real.Project("h_pim_real_tht", "piminus_tht", cut+"piminus_HitTag==1");
	double pim_real = h_pim_real_tht->GetEntries();

	TH1D * h_pip_real_tht = new TH1D("h_pip_real_tht", "tht", 200,0,3.5);
	PiPlus_real.Project("h_pip_real_tht", "piplus_tht", cut+"piplus_HitTag==1");
	double pip_real = h_pip_real_tht->GetEntries();

	TH1D * h_prot_real_tht = new TH1D("h_prot_real_tht", "tht", 200,0,3.5);
	Proton_real.Project("h_prot_real_tht", "proton_tht", cut+"proton_HitTag==1");
	double prot_real = h_prot_real_tht->GetEntries();

	TH1D * h_aprot_real_tht = new TH1D("h_aprot_real_tht", "tht", 200,0,3.5);
	AntiProton_real.Project("h_aprot_real_tht", "AntiProton_tht", cut+"AntiProton_HitTag==1");
	double aprot_real = h_aprot_real_tht->GetEntries();

	TH1D * h_km_real_tht = new TH1D("h_km_real_tht", "tht", 200,0,3.5);
	KaonMinus_real.Project("h_km_real_tht", "kaonminus_tht", cut+"kaonminus_HitTag==1");
	double km_real = h_km_real_tht->GetEntries();

	TH1D * h_kp_real_tht = new TH1D("h_kp_real_tht", "tht", 200,0,3.5);
	KaonPlus_real.Project("h_kp_real_tht", "kaonplus_tht", cut+"kaonplus_HitTag==1");
	double kp_real = h_kp_real_tht->GetEntries();

	cout << "Particle" << "|" << " reco eff (ideal PR) [%]" << "|" << "reco eff (real PR) [%]"   << "|" << "difference [%]" << endl;
	cout << "PiMinus"  << "|" << pim_ideal/nevents_ideal*100 << "|" << pim_real/nevents_real*100  << "|" << (1-pim_real/pim_ideal)*100 << endl;
	cout << "PiPlus"  << "|" << pip_ideal/nevents_ideal*100 << "|" << pip_real/nevents_real*100    << "|" << (1-pip_real/pip_ideal)*100 << endl;
	cout << "Proton"  << "|" << prot_ideal/nevents_ideal*100 << "|" << prot_real/nevents_real*100    << "|" << (1-prot_real/prot_ideal)*100 << endl;
	cout << "AntiProton"  << "|" << aprot_ideal/nevents_ideal*100  << "|" << aprot_real/nevents_real*100  << "|" << (1-aprot_real/aprot_ideal)*100 << endl;
	cout << "KMinus"  << "|" << km_ideal/nevents_ideal*100  << "|" <<  km_real/nevents_real*100  << "|" << (1-km_real/km_ideal)*100 << endl;
	cout << "KPlus"  << "|" << kp_ideal/nevents_ideal*100  << "|" << kp_real/nevents_real*100  << "|" << (1-kp_real/kp_ideal)*100 << endl;

	TCanvas * c = new TCanvas("c", "c", 0,0,1600,1000);
	gStyle->SetOptStat(0);


	//plot for ideal
	TH1D * h_reco_eff_ideal = new TH1D("h_reco_eff_ideal", "; Particle Type; Efficiency", 6,0,6);

	h_reco_eff_ideal->Fill("#pi^{-}", pim_ideal/nevents_ideal);
	h_reco_eff_ideal->Fill("#pi^{+}", pip_ideal/nevents_ideal);
	h_reco_eff_ideal->Fill("p", prot_ideal/nevents_ideal);
	h_reco_eff_ideal->Fill("#bar{p}", aprot_ideal/nevents_ideal);
	h_reco_eff_ideal->Fill("K^{-}", km_ideal/nevents_ideal);
	h_reco_eff_ideal->Fill("K^{+}", kp_ideal/nevents_ideal);

	h_reco_eff_ideal->GetYaxis()->SetRangeUser(0,1);
	h_reco_eff_ideal->GetYaxis()->CenterTitle();
	h_reco_eff_ideal->SetLineColor(kBlue);
	h_reco_eff_ideal->GetXaxis()->SetLabelSize(0.06);

	h_reco_eff_ideal->Draw();

	//plot for real
	TH1D * h_reco_eff_real = new TH1D("h_reco_eff_real", "; Particle Type; Efficiency", 6,0,6);

	h_reco_eff_real->Fill("#pi^{-}", pim_real/nevents_real);
	h_reco_eff_real->Fill("#pi^{+}", pip_real/nevents_real);
	h_reco_eff_real->Fill("p", prot_real/nevents_real);
	h_reco_eff_real->Fill("#bar{p}", aprot_real/nevents_real);
	h_reco_eff_real->Fill("K^{-}", km_real/nevents_real);
	h_reco_eff_real->Fill("K^{+}", kp_real/nevents_real);

	h_reco_eff_real->GetYaxis()->SetRangeUser(0,1);
	h_reco_eff_real->SetLineColor(kRed);

	h_reco_eff_real->Draw("SAME");

  	TLegend * legend = new TLegend(0.76,0.80,0.996,0.92, "");
  	legend->AddEntry(h_reco_eff_ideal, "ideal PR", "l");
  	legend->AddEntry(h_reco_eff_real, "real PR", "l");
  	legend->Draw();


  	c->Print("comparison_PR.pdf");
	c->Print("comparison_PR.png");
	c->Print("comparison_PR.root");



}
