

void comparison_day1_and_withgems(TString input_ideal = "", TString input_real=""){
	if(input_ideal=="" || input_real==""){
		cout << "Please insert two file names" << endl;
	}

	//Open input files and get data
	TFile * data_ideal = new TFile(input_ideal, "READ");
	TFile * data_real = new TFile(input_real, "READ");

	TTree * MC_ideal = (TTree*) data_ideal->Get("ntpMC");
	TTree * piminus_ideal = (TTree*) data_ideal->Get("ntpPiMinus");
	TTree * piplus_ideal = (TTree*) data_ideal->Get("ntpPiPlus");
	TTree * proton_ideal = (TTree*) data_ideal->Get("ntpProton");
	TTree * antiProton_ideal = (TTree*) data_ideal->Get("ntpAntiProton");
	TTree * kaonminus_ideal = (TTree*) data_ideal->Get("ntpKaonMinus");
	TTree * kaonplus_ideal = (TTree*) data_ideal->Get("ntpKaonPlus");

	TTree * MC_real = (TTree*) data_real->Get("ntpMC");
	TTree * piminus_real = (TTree*) data_real->Get("ntpPiMinus");
	TTree * piplus_real = (TTree*) data_real->Get("ntpPiPlus");
	TTree * proton_real = (TTree*) data_real->Get("ntpProton");
	TTree * antiProton_real = (TTree*) data_real->Get("ntpAntiProton");
	TTree * kaonminus_real = (TTree*) data_real->Get("ntpKaonMinus");
	TTree * kaonplus_real = (TTree*) data_real->Get("ntpKaonPlus");

	double nevents_ideal = MC_ideal->GetEntriesFast();
	double nevents_real = MC_real->GetEntriesFast();



	TString cut = "McTruthMatch==1 &";
	TString cut_real = "& McTruthMatch==1";


	//Histograms
	// ideal
	TH1D * h_pim_ideal_tht = new TH1D("h_pim_ideal_tht", "tht", 200,0,3.5);
	piminus_ideal->Project("h_pim_ideal_tht", "piminus_tht", cut+"piminus_HitTag==1 & Mother==3122 ");
	int pim_ideal = h_pim_ideal_tht->GetEntries();

	TH1D * h_pip_ideal_tht = new TH1D("h_pip_ideal_tht", "tht", 200,0,3.5);
	piplus_ideal->Project("h_pip_ideal_tht", "piplus_tht", cut+"piplus_HitTag==1 & (Mother==-3122 || Mother==-3312)");
	int pip_ideal = h_pip_ideal_tht->GetEntries();


	TH1D * h_prot_ideal_tht = new TH1D("h_prot_ideal_tht", "tht", 200,0,3.5);
	proton_ideal->Project("h_prot_ideal_tht", "proton_tht", "Mother==3122 & proton_HitTag==1 & McTruthMatch==1");
	int prot_ideal = h_prot_ideal_tht->GetEntries();

	TH1D * h_aprot_ideal_tht = new TH1D("h_aprot_ideal_tht", "tht", 200,0,3.5);
	antiProton_ideal->Project("h_aprot_ideal_tht", "AntiProton_tht", cut+"AntiProton_HitTag==1 & Mother==-3122");
	int aprot_ideal = h_aprot_ideal_tht->GetEntries();

	TH1D * h_km_ideal_tht = new TH1D("h_km_ideal_tht", "tht", 200,0,3.5);
	kaonminus_ideal->Project("h_km_ideal_tht", "kaonminus_tht", cut+"kaonminus_HitTag==1 & Mother==23314");
	int km_ideal = h_km_ideal_tht->GetEntries();


	// real
	TH1D * h_pim_real_tht = new TH1D("h_pim_real_tht", "tht", 200,0,3.5);
	piminus_real->Project("h_pim_real_tht", "piminus_tht", cut+"piminus_HitTag==1 & Mother==3122 ");
	int pim_real = h_pim_real_tht->GetEntries();

	TH1D * h_pip_real_tht = new TH1D("h_pip_real_tht", "tht", 200,0,3.5);
	piplus_real->Project("h_pip_real_tht", "piplus_tht", cut+"piplus_HitTag==1 & (Mother==-3122 || Mother==-3312)");
	int pip_real = h_pip_real_tht->GetEntries();


	TH1D * h_prot_real_tht = new TH1D("h_prot_real_tht", "tht", 200,0,3.5);
	proton_real->Project("h_prot_real_tht", "proton_tht", "Mother==3122 & proton_HitTag==1 & McTruthMatch==1");
	int prot_real = h_prot_real_tht->GetEntries();

	TH1D * h_aprot_real_tht = new TH1D("h_aprot_real_tht", "tht", 200,0,3.5);
	antiProton_real->Project("h_aprot_real_tht", "AntiProton_tht", cut+"AntiProton_HitTag==1 & Mother==-3122");
	int aprot_real = h_aprot_real_tht->GetEntries();

	TH1D * h_km_real_tht = new TH1D("h_km_real_tht", "tht", 200,0,3.5);
	kaonminus_real->Project("h_km_real_tht", "kaonminus_tht", cut+"kaonminus_HitTag==1 & Mother==23314");
	int km_real = h_km_real_tht->GetEntries();

	cout << "Particle" << "|" << " reco eff (ideal) [%]" << "|" << "reco eff (real) [%]"  << "|" << "difference [%]" << endl;
	cout << "PiMinus"  << "|" << pim_ideal/nevents_ideal*100 << "|" << pim_real/nevents_real*100 << "|" << (1-pim_real/pim_ideal)*100 << endl;
	cout << "PiPlus"  << "|" << pip_ideal/nevents_ideal/2*100 << "|" << pip_real/nevents_real/2*100  << "|" << (1-pip_real/pip_ideal)*100 << endl;
	cout << "Proton"  << "|" << prot_ideal/nevents_ideal*100 << "|" << prot_real/nevents_real*100  << "|" << (1-prot_real/prot_ideal)*100 << endl;
	cout << "Antiproton"  << "|" << aprot_ideal/nevents_ideal*100 << "|" << aprot_real/nevents_real*100 << "|" << (1-aprot_real/aprot_ideal)*100 << endl;
	cout << "KMinus"  << "|" << km_ideal/nevents_ideal*100 << "|" <<  km_real/nevents_real*100  << "|" << (1-km_real/km_ideal)*100 << endl;

	TCanvas * c = new TCanvas("c", "c", 0,0,800,500);
	gStyle->SetOptStat(0);


	//plot for ideal
	TH1D * h_reco_eff_ideal = new TH1D("h_reco_eff_ideal", "; particle type; efficiency", 5,0,5);

	h_reco_eff_ideal->Fill("#pi^{-}", pim_ideal/nevents_ideal);
	h_reco_eff_ideal->Fill("#pi^{+}", pip_ideal/nevents_ideal/2);
	h_reco_eff_ideal->Fill("K^{-}", km_ideal/nevents_ideal);
	h_reco_eff_ideal->Fill("p", prot_ideal/nevents_ideal);
	h_reco_eff_ideal->Fill("#bar{p}", aprot_ideal/nevents_ideal);


	h_reco_eff_ideal->GetYaxis()->SetRangeUser(0,1);
	h_reco_eff_ideal->GetYaxis()->CenterTitle();

	h_reco_eff_ideal->Draw();

	//plot for real
	TH1D * h_reco_eff_real = new TH1D("h_reco_eff_real", "; particle type; efficiency", 5,0,5);

	h_reco_eff_real->Fill("#pi^{-}", pim_real/nevents_real);
	h_reco_eff_real->Fill("#pi^{+}", pip_real/nevents_real/2);
	h_reco_eff_real->Fill("K^{-}", km_real/nevents_real);
	h_reco_eff_real->Fill("p", prot_real/nevents_real);
	h_reco_eff_real->Fill("#bar{p}", aprot_real/nevents_real);


	h_reco_eff_real->GetYaxis()->SetRangeUser(0,1);
	h_reco_eff_real->SetLineColor(kRed);

	h_reco_eff_real->Draw("SAME");

  	TLegend * legend = new TLegend(0.65,0.77,0.889,0.89, "");
  	legend->AddEntry(h_reco_eff_ideal, "day1+GEM", "l");
  	legend->AddEntry(h_reco_eff_real, "day1", "l");
  	legend->Draw();

  	c->Print("comparison_new.pdf");
	c->Print("comparison_new.png");
	c->Print("comparison.root");



}
