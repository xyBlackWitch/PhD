
void pion_kaon_mismatch(TString input=""){


	//*** input file

	TFile * inFile = new TFile(input, "READ");

	//*** Get data from Tree
	TTree * ntpPiMinus = (TTree*) inFile->Get("ntpPiMinus");
	TTree * ntpPiPlus = (TTree*)  inFile->Get("ntpPiPlus");
	TTree * ntpKaonMinus = (TTree*)  inFile->Get("ntpKaonMinus");
	TTree * ntpKaonPlus = (TTree*)  inFile->Get("ntpKaonPlus");


	//*** create histograms to get particle numbers
	TH1D * h_pim = new TH1D("h_pim", "theta distribution; #theta [rad]; counts", 100,0,20);
	TH1D * h_pip = new TH1D("h_pip", "theta distribution; #theta [rad]; counts", 100,0,20);
	TH1D * h_kminus = new TH1D("h_kminus", "theta distribution; #theta [rad]; counts", 100,0,20);
	TH1D * h_kplus = new TH1D("h_kplus", "theta distribution; #theta [rad]; counts", 100,0,20);



	//*** piminus
	ntpPiMinus->Project("h_pim", "piminus_tht", "piminus_HitTag==1");
	float nall_pim = h_pim->GetEntries();

	ntpPiMinus->Project("h_pim", "piminus_tht", "piminus_HitTag==1 & piminus_MC_pdg==-321");
	float nkm_pim = h_pim->GetEntries();

	ntpPiMinus->Project("h_pim", "piminus_tht", "piminus_HitTag==1 & piminus_MC_pdg==321");
	float nkp_pim = h_pim->GetEntries();

	float rate_pim = (nkm_pim + nkp_pim)/nall_pim;

	cout << "#pi^{-}: N(all)" << nall_pim << " N(kminus): " << nkm_pim << " N(kplus): " << nkp_pim << " ratio pion kaon mismatch [%]: " << 100*rate_pim << endl;


	//*** piplus
	ntpPiPlus->Project("h_pip", "piplus_tht", "piplus_HitTag==1");
	float nall_pip = h_pip->GetEntries();

	ntpPiPlus->Project("h_pip", "piplus_tht", "piplus_HitTag==1 & piplus_MC_pdg==-321");
	float nkm_pip = h_pip->GetEntries();

	ntpPiPlus->Project("h_pip", "piplus_tht", "piplus_HitTag==1 & piplus_MC_pdg==321");
	float nkp_pip = h_pip->GetEntries();

	float rate_pip = (nkm_pip + nkp_pip)/nall_pip;

	cout << "#pi^{+}: N(all)" << nall_pip << " N(kminus): " << nkm_pip << " N(kplus): " << nkp_pip << " ratio pion kaon mismatch [%]: " << 100*rate_pip << endl;



	//*** kaonminus
	ntpKaonMinus->Project("h_kminus", "kaonminus_tht", "kaonminus_HitTag==1");
	float nall_kminus = h_kminus->GetEntries();

	ntpKaonMinus->Project("h_kminus", "kaonminus_tht", "kaonminus_HitTag==1 & kaonminus_MC_pdg==-211");
	float nkm_kminus = h_kminus->GetEntries();


	ntpKaonMinus->Project("h_kminus", "kaonminus_tht", "kaonminus_HitTag==1 & kaonminus_MC_pdg==211");
	float nkp_kminus = h_kminus->GetEntries();

	float rate_kminus = (nkm_kminus + nkp_kminus)/nall_kminus;

	cout << "K^{-}: N(all)" << nall_kminus << " N(piminus): " << nkm_kminus << " N(piplus): " << nkp_kminus << " ratio pion kaon mismatch [%]: " << 100*rate_kminus << endl;



	//*** kaonplus
	ntpKaonPlus->Project("h_kplus", "kaonplus_tht", "kaonplus_HitTag==1");
	float nall_kplus = h_kplus->GetEntries();

	ntpKaonPlus->Project("h_kplus", "kaonplus_tht", "kaonplus_HitTag==1 & kaonplus_MC_pdg==-211");
	float nkm_kplus = h_kplus->GetEntries();

	ntpKaonPlus->Project("h_kplus", "kaonplus_tht", "kaonplus_HitTag==1 & kaonplus_MC_pdg==211");
	float nkp_kplus = h_kplus->GetEntries();

	float rate_kplus = (nkm_kplus + nkp_kplus)/nall_kplus;

	cout << "K^{+}: N(all)" << nall_kplus << " N(piplus): " << nkm_kplus << " N(piplus): " << nkp_kplus << " ratio pion kaon mismatch [%]: " << 100*rate_kplus << endl;


}
