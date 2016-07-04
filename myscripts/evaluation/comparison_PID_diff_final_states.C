#include <fstream>
#include <iostream>
#include "/home/ikp1/puetz/panda/PandaSoftware/pandaroot/trunk/source/macro/setPandaStyle.C"


void mycolor(int colornumber, TH1* &hist){

	if(colornumber==0) {hist->SetLineColor(kBlack);}
	else if (colornumber==1) {hist->SetLineColor(kBlue);}
	else if (colornumber==2) {hist->SetLineColor(kRed);}
	else if (colornumber==3) {hist->SetLineColor(kGreen);}
	else if (colornumber==4) {hist->SetLineColor(kViolet);}
	else cout << "no color match!" << endl;

}

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

void comparison_PID_diff_final_states(TString input_list = ""){
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
	TLegend * legend = new TLegend(0.76,0.80,0.996,0.92, "");

	TH1D* h_piminus = new TH1D("h_piminus", "reconstruction efficiencies; PID selection; Efficiency", 6,0,6);
	TH1D* h_piplus = new TH1D("h_piplus", "reconstruction efficiencies; PID selection; Efficiency", 6,0,6);
	TH1D* h_proton = new TH1D("h_proton", "reconstruction efficiencies; PID selection; Efficiency", 6,0,6);
	TH1D* h_antiproton = new TH1D("h_antiproton", "reconstruction efficiencies; PID selection; Efficiency", 6,0,6);
	TH1D* h_kminus = new TH1D("h_kminus", "reconstruction efficiencies; PID selection; Efficiency", 6,0,6);


	while (!filelist.eof()){

		filelist.getline(line, sizeof(line));
		if (strcmp(line, "") == 0) continue;

		TString pid(line);
		unwrap(pid,true);

		TFile * file = new TFile(line, "READ");

		//***Get Data
		TTree * ntpMC = (TTree*) file->Get("ntpMC");
		TTree * ntpPiMinus = (TTree*) file->Get("ntpPiMinus");
		TTree * ntpPiPlus = (TTree*) file->Get("ntpPiPlus");
		TTree * ntpProton = (TTree*) file->Get("ntpProton");
		TTree * ntpAntiProton = (TTree*) file->Get("ntpAntiProton");
		TTree * ntpKaonMinus = (TTree*) file->Get("ntpKaonMinus");


		double nevents = ntpMC->GetEntriesFast();

		//***Histograms

		TH1D * h_pim = new TH1D("h_pim", "tht", 200,0,10);
		ntpPiMinus->Project("h_pim", "piminus_tht", "Mother==3122 & McTruthMatch & piminus_HitTag==1");
		double pim = h_pim->GetEntries();

		TH1D * h_pip_tht = new TH1D("h_pip_tht", "tht", 200,0,10);
		ntpPiPlus->Project("h_pip_tht", "piplus_tht", "Mother==-3122 & McTruthMatch & piplus_HitTag==1");
		double pip = h_pip_tht->GetEntries();

		TH1D * h_pip2_tht = new TH1D("h_pip2_tht", "tht", 200,0,10);
		ntpPiPlus->Project("h_pip2_tht", "piplus_tht", "Mother==-3312 & McTruthMatch & piplus_HitTag==1");
		double pip2 = h_pip2_tht->GetEntries();

		TH1D * h_prot_tht = new TH1D("h_prot_tht", "tht", 200,0,10);
		ntpProton->Project("h_prot_tht", "proton_tht", "MC_Mother_PDG==3122 & McTruthMatch==1 & proton_HitTag==1");
		double prot = h_prot_tht->GetEntries();

		TH1D * h_aprot_tht = new TH1D("h_aprot_tht", "tht", 200,0,10);
		ntpAntiProton->Project("h_aprot_tht", "AntiProton_tht", "Mother==-3122 & McTruthMatch & AntiProton_HitTag==1");
		double aprot = h_aprot_tht->GetEntries();

		TH1D * h_km_tht = new TH1D("h_km_tht", "tht", 200,0,10);
		ntpKaonMinus->Project("h_km_tht", "kaonminus_tht", "Mother==23314 & McTruthMatch & kaonminus_HitTag==1");
		double km = h_km_tht->GetEntries();




		//***Fill histograms

		h_piminus->Fill(pid, pim/nevents);
		h_piplus->Fill(pid, (pip+pip2)/2/nevents);
		h_proton->Fill(pid, prot/nevents);
		h_antiproton->Fill(pid, aprot/nevents);
		h_kminus->Fill(pid, km/nevents);

		cout << "PID selection: " << pid << endl;
		cout << "Particle" << "|" << " reco eff (ideal) [%]" << endl;
		cout << "PiMinus"  << "|" << pim/nevents*100 << endl;
		cout << "PiPlus"  << "|" << (pip+pip2)/2/nevents*100 << endl;
		cout << "Proton"  << "|" << prot/nevents*100 << endl;
		cout << "Antiproton"  << "|" << aprot/nevents*100 << endl;
		cout << "KMinus"  << "|" << km/nevents*100 << endl;

		cout << " " << endl;

	}

	h_piminus->GetYaxis()->SetRangeUser(0,1);
	mycolor(0, h_piminus);
	mycolor(1, h_piplus);
	mycolor(2, h_proton);
	mycolor(3, h_antiproton);
	mycolor(4, h_kminus);

  	legend->AddEntry(h_piminus, "#pi^{-}", "l");
	legend->AddEntry(h_piplus, "#pi^{+}", "l");
	legend->AddEntry(h_proton, "p", "l");
	legend->AddEntry(h_antiproton, "#bar{p}", "l");
	legend->AddEntry(h_kminus, "K^{-}", "l");


	h_piminus->Draw("PL");
	h_piplus->Draw("SAME PL");
	h_proton->Draw("SAME PL");
	h_antiproton->Draw("SAME PL");
	h_kminus->Draw("SAME PL");

	legend->Draw();



  	c->Print("comparison_PID_FS.pdf");
	c->Print("comparison_PID_FS.png");
	c->Print("comparison_PID_FS.root");

}
