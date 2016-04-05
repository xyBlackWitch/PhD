

void comparison_GenFit(TString input_GenFit1 = "", TString input_GenFit2=""){
	if(input_GenFit1=="" || input_GenFit2==""){
		cout << "Please insert two file names" << endl;
	}

	//Open input files and get data
	TFile * data_GenFit1 = new TFile(input_GenFit1, "READ");
	TFile * data_GenFit2 = new TFile(input_GenFit2, "READ");

	TTree * MC_GenFit1 = (TTree*) data_GenFit1->Get("ntpMC");
	TTree * piminus_GenFit1 = (TTree*) data_GenFit1->Get("ntpPiMinus");
	TTree * piplus_GenFit1 = (TTree*) data_GenFit1->Get("ntpPiPlus");
	TTree * proton_GenFit1 = (TTree*) data_GenFit1->Get("ntpProton");
	TTree * antiProton_GenFit1 = (TTree*) data_GenFit1->Get("ntpAntiProton");
	TTree * kaonminus_GenFit1 = (TTree*) data_GenFit1->Get("ntpKaonMinus");

	TTree * MC_GenFit2 = (TTree*) data_GenFit2->Get("ntpMC");
	TTree * piminus_GenFit2 = (TTree*) data_GenFit2->Get("ntpPiMinus");
	TTree * piplus_GenFit2 = (TTree*) data_GenFit2->Get("ntpPiPlus");
	TTree * proton_GenFit2 = (TTree*) data_GenFit2->Get("ntpProton");
	TTree * antiProton_GenFit2 = (TTree*) data_GenFit2->Get("ntpAntiProton");
	TTree * kaonminus_GenFit2 = (TTree*) data_GenFit2->Get("ntpKaonMinus");

	double nevents_GenFit1 = MC_GenFit1->GetEntriesFast();
	double nevents_GenFit2 = MC_GenFit2->GetEntriesFast();

	cout << nevents_GenFit2 << endl;

	TString cut = "McTruthMatch==1 &";


	//Histograms
	// GenFit1
	TH1D * h_pim_GenFit1_tht = new TH1D("h_pim_GenFit1_tht", "tht", 200,0,3.5);
	piminus_GenFit1->Project("h_pim_GenFit1_tht", "piminus_tht", cut+"piminus_HitTag==1 & Mother==3122");
	int pim_GenFit1 = h_pim_GenFit1_tht->GetEntries();

	TH1D * h_pip_GenFit1_tht = new TH1D("h_pip_GenFit1_tht", "tht", 200,0,3.5);
	piplus_GenFit1->Project("h_pip_GenFit1_tht", "piplus_tht", cut+"piplus_HitTag==1 & Mother==-3122");
	int pip_GenFit1 = h_pip_GenFit1_tht->GetEntries();

	TH1D * h_pip_GenFit1_tht2 = new TH1D("h_pip_GenFit1_tht2", "tht", 200,0,3.5);
	piplus_GenFit1->Project("h_pip_GenFit1_tht2", "piplus_tht", cut+"piplus_HitTag==1 & Mother==-3312");
	int pip2_GenFit1 = h_pip_GenFit1_tht2->GetEntries();

	TH1D * h_prot_GenFit1_tht = new TH1D("h_prot_GenFit1_tht", "tht", 200,0,3.5);
	proton_GenFit1->Project("h_prot_GenFit1_tht", "proton_tht", "MC_Mother_PDG==3122 & proton_HitTag==1 & McTruthMatch==1");
	int prot_GenFit1 = h_prot_GenFit1_tht->GetEntries();

	TH1D * h_aprot_GenFit1_tht = new TH1D("h_aprot_GenFit1_tht", "tht", 200,0,3.5);
	antiProton_GenFit1->Project("h_aprot_GenFit1_tht", "AntiProton_tht", cut+"AntiProton_HitTag==1 & Mother==-3122");
	int aprot_GenFit1 = h_aprot_GenFit1_tht->GetEntries();

	TH1D * h_km_GenFit1_tht = new TH1D("h_km_GenFit1_tht", "tht", 200,0,3.5);
	kaonminus_GenFit1->Project("h_km_GenFit1_tht", "kaonminus_tht", cut+"kaonminus_HitTag==1 & Mother==23314");
	int km_GenFit1 = h_km_GenFit1_tht->GetEntries();



	//GenFit2
	TH1D * h_pim_GenFit2_tht = new TH1D("h_pim_GenFit2_tht", "tht", 200,0,3.5);
	piminus_GenFit2->Project("h_pim_GenFit2_tht", "piminus_tht", cut+"piminus_HitTag==1 & Mother==3122");
	double pim_GenFit2 = h_pim_GenFit2_tht->GetEntries();

	TH1D * h_pip_GenFit2_tht = new TH1D("h_pip_GenFit2_tht", "tht", 200,0,3.5);
	piplus_GenFit2->Project("h_pip_GenFit2_tht", "piplus_tht", cut+"piplus_HitTag==1 & Mother==-3122");
	double pip_GenFit2 = h_pip_GenFit2_tht->GetEntries();

	TH1D * h_pip_GenFit2_tht2 = new TH1D("h_pip_GenFit2_tht2", "tht", 200,0,3.5);
	piplus_GenFit2->Project("h_pip_GenFit2_tht2", "piplus_tht", cut+"piplus_HitTag==1 & Mother==-3312");
	double pip2_GenFit2 = h_pip_GenFit2_tht2->GetEntries();

	TH1D * h_prot_GenFit2_tht = new TH1D("h_prot_GenFit2_tht", "tht", 200,0,3.5);
	proton_GenFit2->Project("h_prot_GenFit2_tht", "proton_tht", "MC_Mother_PDG==3122 & proton_HitTag==1 & McTruthMatch==1");
	double prot_GenFit2 = h_prot_GenFit2_tht->GetEntries();

	TH1D * h_aprot_GenFit2_tht = new TH1D("h_aprot_GenFit2_tht", "tht", 200,0,3.5);
	antiProton_GenFit2->Project("h_aprot_GenFit2_tht", "AntiProton_tht", cut+"AntiProton_HitTag==1 & Mother==-3122");
	double aprot_GenFit2 = h_aprot_GenFit2_tht->GetEntries();

	TH1D * h_km_GenFit2_tht = new TH1D("h_km_GenFit2_tht", "tht", 200,0,3.5);
	kaonminus_GenFit2->Project("h_km_GenFit2_tht", "kaonminus_tht", cut+"kaonminus_HitTag==1 & Mother==23314");
	double km_GenFit2 = h_km_GenFit2_tht->GetEntries();



	cout << "Particle" << "|" << " reco eff (GenFit1) [%]" << "|" << "reco eff (GenFit2) [%]"  << "|" << "difference [%]" << endl;
	cout << "PiMinus"  << "|" << pim_GenFit1/nevents_GenFit1*100 << "|" << pim_GenFit2/nevents_GenFit2*100 << "|" << (1-pim_GenFit2/pim_GenFit1)*100 << endl;
	cout << "PiPlus(alam)"  << "|" << pip_GenFit1/nevents_GenFit1*100 << "|" << pip_GenFit2/nevents_GenFit2*100  << "|" << (1-pip_GenFit2/pip_GenFit1)*100 << endl;
	cout << "PiPlus(Xi+)"  << "|" << pip2_GenFit1/nevents_GenFit1*100 << "|" << pip2_GenFit2/nevents_GenFit2*100  << "|" << (1-pip2_GenFit2/pip2_GenFit1)*100 << endl;
	cout << "Proton"  << "|" << prot_GenFit1/nevents_GenFit1*100 << "|" << prot_GenFit2/nevents_GenFit2*100  << "|" << (1-prot_GenFit2/prot_GenFit1)*100 << endl;
	cout << "Antiproton"  << "|" << aprot_GenFit1/nevents_GenFit1*100 << "|" << aprot_GenFit2/nevents_GenFit2*100 << "|" << (1-aprot_GenFit2/aprot_GenFit1)*100 << endl;
	cout << "KMinus"  << "|" << km_GenFit1/nevents_GenFit1*100 << "|" <<  km_GenFit2/nevents_GenFit2*100  << "|" << (1-km_GenFit2/km_GenFit1)*100 << endl;

	TCanvas * c = new TCanvas("c", "c", 0,0,800,500);
	gStyle->SetOptStat(0);


	//plot for GenFit1
	TH1D * h_reco_eff_GenFit1 = new TH1D("h_reco_eff_GenFit1", "; particle type; efficiency", 5,0,5);

	h_reco_eff_GenFit1->Fill("#pi^{-}", pim_GenFit1/nevents_GenFit1);
	h_reco_eff_GenFit1->Fill("#pi^{+}", (pip_GenFit1+pip2_GenFit1)/2/nevents_GenFit1);
	h_reco_eff_GenFit1->Fill("p", prot_GenFit1/nevents_GenFit1);
	h_reco_eff_GenFit1->Fill("#bar{p}", aprot_GenFit1/nevents_GenFit1);
	h_reco_eff_GenFit1->Fill("K^{-}", km_GenFit1/nevents_GenFit1);

	h_reco_eff_GenFit1->GetYaxis()->SetRangeUser(0,1);
	h_reco_eff_GenFit1->GetYaxis()->CenterTitle();

	h_reco_eff_GenFit1->Draw();

	//plot for GenFit2
	TH1D * h_reco_eff_GenFit2 = new TH1D("h_reco_eff_GenFit2", "; particle type; efficiency", 5,0,5);

	h_reco_eff_GenFit2->Fill("#pi^{-}", pim_GenFit2/nevents_GenFit2);
	h_reco_eff_GenFit2->Fill("#pi^{+}", (pip_GenFit2+pip2_GenFit2)/2/nevents_GenFit2);
	h_reco_eff_GenFit2->Fill("p", prot_GenFit2/nevents_GenFit2);
	h_reco_eff_GenFit2->Fill("#bar{p}", aprot_GenFit2/nevents_GenFit2);
	h_reco_eff_GenFit2->Fill("K^{-}", km_GenFit2/nevents_GenFit2);

	h_reco_eff_GenFit2->GetYaxis()->SetRangeUser(0,1);
	h_reco_eff_GenFit2->SetLineColor(kRed);

	h_reco_eff_GenFit2->Draw("SAME");

  	TLegend * legend = new TLegend(0.65,0.77,0.889,0.89, "");
  	legend->AddEntry(h_reco_eff_GenFit1, "GenFit1", "l");
  	legend->AddEntry(h_reco_eff_GenFit2, "GenFit2", "l");
  	legend->Draw();




}
