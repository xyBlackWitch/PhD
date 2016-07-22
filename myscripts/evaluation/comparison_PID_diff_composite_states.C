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

void comparison_PID_diff_composite_states(TString input_list = "", TString particle = "lambda0"){
	if(input_list==""){
		cout << "Please insert txt-file name" << endl;
	}

	std::ifstream filelist(input_list);

	if (!filelist.is_open()){
		cerr << "failed in open." << endl;
		exit();
	}

	char line[10000];
	TString ntp;
	TString name;

	setPandaStyle();
	TCanvas * c = new TCanvas("c", "c", 0,0,1600,1000);

	TString title = TString::Format("reconstruction efficiencies for %s; PID selection; Efficiency", particle.Data());

	TH1D* h_composite = new TH1D("h_composite", title, 6,0,6);


	//*** Which particle is chosen?
	if(particle=="lambda"){
		ntp="ntpLambda0";
		name="Lambda0";
	}
	else if (particle=="antilambda"){
		ntp="ntpAntiLambda0";
		name="antiLambda0";
	}
	else if (particle=="xiplus"){
		ntp="ntpXiPlus";
		name="xiplus";
	}
	else if (particle=="ximinus1820"){
		ntp="ntpXiMinus1820";
		name="XiMinus";
	}
	else{
		cout << "Please insert a valid particle! (lambda, antilambda, xiplus or ximinus1820)" << endl;
	}

	cout << "PID selection" << "|" << " reco eff (ideal) [%]" << endl;

	while (!filelist.eof()){

		filelist.getline(line, sizeof(line));
		if (strcmp(line, "") == 0) continue;

		TString pid(line);
		unwrap(pid,true);

		TFile * file = new TFile(line, "READ");

		//***Get Data
		TTree * ntpMC = (TTree*) file->Get("ntpMC");
		TTree * ntpSys = (TTree*) file->Get(ntp);

		double nevents = ntpMC->GetEntriesFast();

		//***Histograms

		TString cut = TString::Format("VtxFit_HowGood==1 & MassFit_prob>0.01 & %s_HitTag==1 & McTruthMatch", name.Data());

		TH1D * h_cand = new TH1D("h_cand", "tht", 200,0,10);
		ntpSys->Project("h_cand", name+"_tht", cut);
		double sys = h_cand->GetEntries();


		//***Fill histograms

		h_composite->Fill(pid, sys/nevents);

		cout << pid  << " |" << sys/nevents*100 << endl;

	}

	h_composite->GetYaxis()->SetRangeUser(0,1);
	h_composite->Draw("PL");

	c->Print("comparison_PID_"+particle+".pdf");
	c->Print("comparison_PID_"+particle+".png");
	c->Print("comparison_PID_"+particle+".root");

}
