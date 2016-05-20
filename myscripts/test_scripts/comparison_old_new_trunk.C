

void comparison_old_new_trunk(TString input_old = "", TString input_new=""){
	if(input_old=="" || input_new==""){
		cout << "Please insert two file names" << endl;
	}

	//Open input files and get data
	TFile * data_old = new TFile(input_old, "READ");
	TFile * data_new = new TFile(input_new, "READ");

	TTree * MC_old = (TTree*) data_old->Get("ntpMC");
	TTree * piminus_old = (TTree*) data_old->Get("ntpPiMinus");
	TTree * piplus_old = (TTree*) data_old->Get("ntpPiPlus");
	TTree * proton_old = (TTree*) data_old->Get("ntpProton");
	TTree * antiProton_old = (TTree*) data_old->Get("ntpAntiProton");
	TTree * kaonminus_old = (TTree*) data_old->Get("ntpKaonMinus");


	TTree * MC_new = (TTree*) data_new->Get("ntpMC");
	TTree * piminus_new = (TTree*) data_new->Get("ntpPiMinus");
	TTree * piplus_new = (TTree*) data_new->Get("ntpPiPlus");
	TTree * proton_new = (TTree*) data_new->Get("ntpProton");
	TTree * antiProton_new = (TTree*) data_new->Get("ntpAntiProton");
	TTree * kaonminus_new = (TTree*) data_new->Get("ntpKaonMinus");


	double nevents_old = MC_old->GetEntriesFast();
	double nevents_new = MC_new->GetEntriesFast();

	cout << "Nevents old: " << nevents_old << " and Nevents new: " << nevents_new << endl;

	TString cut = "McTruthMatch==1 & ";
	TString cut_new = "& McTruthMatch==1";


	//Histograms
	// old
	TH1D * h_pim_old_tht = new TH1D("h_pim_old_tht", "tht", 200,0,3.5);
	piminus_old->Project("h_pim_old_tht", "piminus_tht", cut+"piminus_HitTag==1 & Mother==3122");
	int pim_old = h_pim_old_tht->GetEntries();

	TH1D * h_pip_old_tht = new TH1D("h_pip_old_tht", "tht", 200,0,3.5);
	piplus_old->Project("h_pip_old_tht", "piplus_tht", cut+"piplus_HitTag==1 & Mother==-3122");
	int pip_old = h_pip_old_tht->GetEntries();

	TH1D * h_pip_old_tht2 = new TH1D("h_pip_old_tht2", "tht", 200,0,3.5);
	piplus_old->Project("h_pip_old_tht2", "piplus_tht", cut+"piplus_HitTag==1 & Mother==-3312");
	int pip_old2 = h_pip_old_tht2->GetEntries();

	TH1D * h_prot_old_tht = new TH1D("h_prot_old_tht", "tht", 200,0,3.5);
	proton_old->Project("h_prot_old_tht", "proton_tht", "MC_Mother_PDG==3122 & proton_HitTag==1 & McTruthMatch==1");
	int prot_old = h_prot_old_tht->GetEntries();

	TH1D * h_aprot_old_tht = new TH1D("h_aprot_old_tht", "tht", 200,0,3.5);
	antiProton_old->Project("h_aprot_old_tht", "AntiProton_tht", cut+"AntiProton_HitTag==1 & Mother==-3122");
	int aprot_old = h_aprot_old_tht->GetEntries();

	TH1D * h_km_old_tht = new TH1D("h_km_old_tht", "tht", 200,0,3.5);
	kaonminus_old->Project("h_km_old_tht", "kaonminus_tht", cut+"kaonminus_HitTag==1 & Mother==23314");
	int km_old = h_km_old_tht->GetEntries();



	//new
	TH1D * h_pim_new_tht = new TH1D("h_pim_new_tht", "tht", 200,0,3.5);
	piminus_new->Project("h_pim_new_tht", "piminus_tht", cut+"piminus_HitTag==1 & Mother==3122");
	double pim_new = h_pim_new_tht->GetEntries();

	TH1D * h_pip_new_tht = new TH1D("h_pip_new_tht", "tht", 200,0,3.5);
	piplus_new->Project("h_pip_new_tht", "piplus_tht", cut+"piplus_HitTag==1 & Mother==-3312");
	double pip_new = h_pip_new_tht->GetEntries();

	TH1D * h_pip_new_tht2 = new TH1D("h_pip_new_tht2", "tht", 200,0,3.5);
	piplus_new->Project("h_pip_new_tht2", "piplus_tht", cut+"piplus_HitTag==1 & Mother==-3122");
	double pip_new2 = h_pip_new_tht2->GetEntries();

	TH1D * h_prot_new_tht = new TH1D("h_prot_new_tht", "tht", 200,0,3.5);
	proton_new->Project("h_prot_new_tht", "proton_tht", "MC_Mother_PDG==3122 & McTruthMatch==1 & proton_HitTag==1");
	double prot_new = h_prot_new_tht->GetEntries();

	TH1D * h_aprot_new_tht = new TH1D("h_aprot_new_tht", "tht", 200,0,3.5);
	antiProton_new->Project("h_aprot_new_tht", "AntiProton_tht", cut+"AntiProton_HitTag==1 & Mother==-3122");
	double aprot_new = h_aprot_new_tht->GetEntries();

	TH1D * h_km_new_tht = new TH1D("h_km_new_tht", "tht", 200,0,3.5);
	kaonminus_new->Project("h_km_new_tht", "kaonminus_tht", cut+"kaonminus_HitTag==1 & Mother==23314");
	double km_new = h_km_new_tht->GetEntries();


	cout << "Particle" << "|" << " reco eff (old) [%]" << "|" << "reco eff (new) [%]"  << "|" << "difference [%]" << endl;
	cout << "PiMinus"  << "|" << pim_old/nevents_old*100 << "|" << pim_new/nevents_new*100 << "|" << (1-pim_new/pim_old)*100 << endl;
	cout << "PiPlus"  << "|" << (pip_old+pip_old2)/nevents_old*50 << "|" << (pip_new+pip_new2)/nevents_new*50  << "|" << (1-(pip_new+pip_new2)/(pip_old+pip_old2))*100 << endl;
	cout << "Proton"  << "|" << prot_old/nevents_old*100 << "|" << prot_new/nevents_new*100  << "|" << (1-prot_new/prot_old)*100 << endl;
	cout << "Antiproton"  << "|" << aprot_old/nevents_old*100 << "|" << aprot_new/nevents_new*100 << "|" << (1-aprot_new/aprot_old)*100 << endl;
	cout << "KMinus"  << "|" << km_old/nevents_old*100 << "|" <<  km_new/nevents_new*100  << "|" << (1-km_new/km_old)*100 << endl;

	TCanvas * c = new TCanvas("c", "c", 0,0,800,500);
	gStyle->SetOptStat(0);


	//plot for old
	TH1D * h_reco_eff_old = new TH1D("h_reco_eff_old", "; particle type; efficiency", 5,0,5);

	h_reco_eff_old->Fill("#pi^{-}", pim_old/nevents_old);
	h_reco_eff_old->Fill("#pi^{+}", (pip_old+pip_old2)/2/nevents_old);
	h_reco_eff_old->Fill("p", prot_old/nevents_old);
	h_reco_eff_old->Fill("#bar{p}", aprot_old/nevents_old);
	h_reco_eff_old->Fill("K^{-}", km_old/nevents_old);

	h_reco_eff_old->GetYaxis()->SetRangeUser(0,1);
	h_reco_eff_old->GetYaxis()->CenterTitle();

	h_reco_eff_old->Draw();

	//plot for new
	TH1D * h_reco_eff_new = new TH1D("h_reco_eff_new", "; particle type; efficiency", 5,0,5);

	h_reco_eff_new->Fill("#pi^{-}", pim_new/nevents_new);
	h_reco_eff_new->Fill("#pi^{+}", (pip_new+pip_new2)/2/nevents_new);
	h_reco_eff_new->Fill("p", prot_new/nevents_new);
	h_reco_eff_new->Fill("#bar{p}", aprot_new/nevents_new);
	h_reco_eff_new->Fill("K^{-}", km_new/nevents_new);

	h_reco_eff_new->GetYaxis()->SetRangeUser(0,1);
	h_reco_eff_new->SetLineColor(kRed);

	h_reco_eff_new->Draw("SAME");

  	TLegend * legend = new TLegend(0.65,0.77,0.889,0.89, "");
  	legend->AddEntry(h_reco_eff_old, "old trunk", "l");
  	legend->AddEntry(h_reco_eff_new, "new trunk", "l");
  	legend->Draw();





}
