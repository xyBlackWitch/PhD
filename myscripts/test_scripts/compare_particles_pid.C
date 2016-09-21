

void compare_particles_pid(TString file_all="", TString file_ideal=""){

	TFile * data_all = new TFile(file_all, "READ");
	TFile * data_ideal = new TFile(file_ideal, "READ");

	TTree * all = (TTree*) data_all->Get("ntpPiMinus");
	TTree * ideal = (TTree*) data_ideal->Get("ntpPiMinus");

	int nevents = 100;

	for(int i=0; i<nevents; i++){

		all->GetEntry(i);

		double evt_all = all->GetLeaf("ev")->GetValue(0);
		double cand_all = all->GetLeaf("cand")->GetValue(0);
		double ppi_all = all->GetLeaf("piminus_p")->GetValue(0);
		double hittag_all = all->GetLeaf("piminus_HitTag")->GetValue(0);

		if (hittag_all) cout << "All : Pion from event: " << evt_all << " Candidate: " << cand_all << " with p= " << ppi_all << endl;

		ideal->GetEntry(i);

		double evt_ideal = ideal->GetLeaf("ev")->GetValue(0);
		double cand_ideal = ideal->GetLeaf("cand")->GetValue(0);
		double ppi_ideal = ideal->GetLeaf("piminus_p")->GetValue(0);

		double hittag_ideal = ideal->GetLeaf("piminus_HitTag")->GetValue(0);

		if (hittag_ideal) cout << "Ideal : Pion from event: " << evt_ideal << " Candidate: " << cand_ideal << " with p= " << ppi_ideal << endl;

	}

}
