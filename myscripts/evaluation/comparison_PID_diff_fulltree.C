#include <fstream>
#include <iostream>
#include "/home/ikp1/puetz/panda/PandaSoftware/pandaroot/trunk/source/macro/setPandaStyle.C"


void unwrap(TString &name, int secondeunderscore=false){

	//**remove the file ending
	int end = name.First(".");
	name.Remove(end);

	//**remove prefix
	int underscore = name.First("_");
	name.Remove(0,underscore+1);

	//** check if file name has second underscore und remove it
	if(secondeunderscore){
		int secunder = name.First("_");
		name.Remove(0,secunder+1);
	}

}

void comparison_PID_diff_fulltree(TString input_list = ""){
	if(input_list==""){
		cout << "Please insert txt-file name" << endl;
	}

	std::ifstream filelist(input_list);

	if (!filelist.is_open()){
		cerr << "failed in open." << endl;
		exit();
	}

	char line[10000];

	setPandaStyle();
	TCanvas * c = new TCanvas("c", "c", 0,0,1600,1000);

	TH1D* h_fulltree = new TH1D("h_fulltree", "reconstruction efficiencies; PID selection; Efficiency [%]", 6,0,6);

	cout << "PID selection" << "|" << " reco eff (ideal) [%]" << endl;

	while (!filelist.eof()){

		filelist.getline(line, sizeof(line));
		if (strcmp(line, "") == 0) continue;

		TString pid(line);
		unwrap(pid,true);

		TFile * file = new TFile(line, "READ");

		//***Get Data
		TTree * ntpMC = (TTree*) file->Get("ntpMC");
		TTree * ntpSys = (TTree*) file->Get("ntpXiSys");

		double nevents = ntpMC->GetEntriesFast();

		//***Histograms

		TH1D * h_sys = new TH1D("h_sys", "tht", 200,0,10);
		ntpSys->Project("h_sys", "XiSys_tht", "4CFit_prob>0.01 & McTruthMatch");
		double sys = h_sys->GetEntries();


		//***Fill histograms

		h_fulltree->Fill(pid, sys/nevents*100);

		cout << pid  << " |" << sys/nevents*100 << endl;

	}

	h_fulltree->GetYaxis()->SetRangeUser(0,3);
	h_fulltree->Draw("PL");

	c->Print("comparison_PID_fulltree.pdf");
	c->Print("comparison_PID_fulltree.png");
	c->Print("comparison_PID_fulltree.root");

}
