
void combine_leaving_GEM_hits(){


	TFile * inPim = new TFile("./pions/piminus_withgem_test_output_ana.root", "READ");
	TFile * inPip = new TFile("./pions/piplus_withgem_test_output_ana.root", "READ");
	TFile * inProt = new TFile("./protons/proton_withgem_test_output_ana.root", "READ");
	TFile * inaProt = new TFile("./antiprotons/antiproton_withgem_test_output_ana.root", "READ");
	TFile * inkaonm = new TFile("./kaons/kaonminus_withgem_test_output_ana.root", "READ");
	TFile * inkaonp = new TFile("./kaons/kaonplus_withgem_test_output_ana.root", "READ");

	TTree * treePim = (TTree*) inPim->Get("ntpPiMinus");
	TTree * treePip = (TTree*) inPip->Get("ntpPiPlus");
	TTree * treeProt = (TTree*) inProt->Get("ntpProton");
	TTree * treeaProt = (TTree*) inaProt->Get("ntpAntiProton");
	TTree * treekaonm = (TTree*) inkaonm->Get("ntpKaonMinus");
	TTree * treekaonp = (TTree*) inkaonp->Get("ntpKaonPlus");

	int nevts=10000;

	TH1D * h = new TH1D("h", "particles per event leaving hit in GEM planes; #particles; counts", 7,0,7);
	h->GetXaxis()->SetLabelSize(0.05);
	h->GetXaxis()->SetTitleSize(0.05);

	h->GetYaxis()->SetLabelSize(0.05);
	h->GetYaxis()->SetTitleSize(0.05);


	for(int i=0; i<nevts; i++){

		int sum=0;
		//***** PiMinus *****
		treePim->GetEntry(i);
		int countpim = treePim->GetLeaf("GemHit")->GetValue(0);

		sum += countpim;


		//***** Piplus *****
		treePip->GetEntry(i);
		int countpip = treePip->GetLeaf("GemHit")->GetValue(0);

		sum += countpip;


		//***** Proton *****
		treeProt->GetEntry(i);
		int countProt = treeProt->GetLeaf("GemHit")->GetValue(0);

		sum += countProt;


		//***** AntiProton *****
		treeaProt->GetEntry(i);
		int countaProt = treeaProt->GetLeaf("GemHit")->GetValue(0);

		sum += countaProt;


		//***** kaonMinus *****
		treekaonm->GetEntry(i);
		int countkaonm = treekaonm->GetLeaf("GemHit")->GetValue(0);

		sum += countkaonm;


		//***** kaonplus *****
		treekaonp->GetEntry(i);
		int countkaonp = treekaonp->GetLeaf("GemHit")->GetValue(0);

		sum += countkaonp;

		h->Fill(sum);

	}
	TCanvas * c = new TCanvas("c", "c", 0,0,1000,600);
	h->Draw();
	c->Print("particles_leaving_GemHits_boxgen.root");
	c->Print("particles_leaving_GemHits_boxgen.png");


}
