#include <fstream>
#include <iostream>
#include "/home/ikp1/puetz/panda/PandaSoftware/pandaroot/trunk/source/macro/setPandaStyle.C"


void mycolor(int colornumber, TH1* &hist){

	hist->SetLineColor(kBlack);

	if(colornumber==0) {
		hist->SetFillColor(kYellow);
		hist->SetFillStyle(1001);
	}
	else if (colornumber==1) {
		hist->SetFillColor(kBlue);
		hist->SetFillStyle(3004);
	}
	else if (colornumber==2) {
		hist->SetFillStyle(3001);
		hist->SetFillColor(kRed);
	}
	else if (colornumber==3) {
		hist->SetFillStyle(1001);
		hist->SetFillColor(kGreen);
	}
	else if (colornumber==4) {
		hist->SetFillStyle(1001);
		hist->SetFillColor(kViolet);
	}
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


void pion_kaon_mismatch(TString input_list = "", TString particle="piminus"){

	//*** check if input list is empty
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
		TString sign;



		//***Create Canvas and Histogram
		setPandaStyle();
		TCanvas * c = new TCanvas("c", "pion kaon mismatch",0,0,1600,1000);
		TLegend * legend = new TLegend(0.68,0.84,0.999,0.96, "");

		TString title = TString::Format("%s; PID selector; counts", particle.Data());
		TH1D * h_all = new TH1D("h_all", title, 5,0,5);
		TH1D * h_kaon = new TH1D("h_kaon", title, 5,0,5);
		TH1D * h_other = new TH1D("h_other", "", 5,0,5);

		//*** Which particle is chosen?
		if(particle=="piminus"){
			ntp="ntpPiMinus";
			sign="-";
		}
		else if (particle=="piplus"){
			ntp="ntpPiPlus";
			sign="";
		}
		else if (particle=="kaonminus"){
			ntp="ntpKaonMinus";
			sign="-";
		}
		else if (particle=="kaonplus"){
			ntp="ntpKaonPlus";
			sign="";
		}
		else{
			cout << "Please insert a valid particle! (piminus, piplus, kaonminus or kaonplus)" << endl;
		}

		while (!filelist.eof()){

			filelist.getline(line, sizeof(line));
			if (strcmp(line, "") == 0) continue;

			TString pid(line);
			unwrap(pid,true);

			if(pid=="idealPID") continue;

			TFile * file = new TFile(line, "READ");

			//** Get Data
			TTree * data = (TTree*) file->Get(ntp);

			if (particle=="piminus" || particle=="piplus"){
				TString cut_all = TString::Format("%s_HitTag==1", particle.Data());
				TString cut_kaon = TString::Format("%s_HitTag==1 & %s_MC_pdg==%s321", particle.Data(), particle.Data(),sign.Data());
				TString cut_other = TString::Format("%s_HitTag==1 & %s_MC_pdg!=%s321 & %s_MC_pdg!=%s211", particle.Data(), particle.Data(),sign.Data(), particle.Data(), sign.Data());
			}
			else{
				TString cut_all = TString::Format("%s_HitTag==1", particle.Data());
				TString cut_kaon = TString::Format("%s_HitTag==1 & %s_MC_pdg==%s211", particle.Data(), particle.Data(),sign.Data());
				TString cut_other = TString::Format("%s_HitTag==1 & %s_MC_pdg!=%s321 & %s_MC_pdg!=%s211", particle.Data(), particle.Data(),sign.Data(), particle.Data(), sign.Data());

			}

			TH1D * h_tht_all = new TH1D("h_tht_all", "", 200,0,10);
			data->Project("h_tht_all", particle+"_tht", cut_all);
			int n_all = h_tht_all->GetEntries();

			TH1D * h_tht_kaon = new TH1D("h_tht_kaon", "", 200,0,10);
			data->Project("h_tht_kaon", particle+"_tht", cut_kaon);
			int n_kaon = h_tht_kaon->GetEntries();

			TH1D * h_tht_other = new TH1D("h_tht_other", "", 200,0,10);
			data->Project("h_tht_other", particle+"_tht", cut_other);
			int n_other = h_tht_other->GetEntries();

			h_all->Fill(pid, n_all);
			h_kaon->Fill(pid,n_kaon);
			h_other->Fill(pid, n_other);



		}

		mycolor(0,h_all);
		mycolor(1,h_kaon);
		mycolor(2,h_other);

		if (particle=="piminus" || particle=="piplus"){
			legend->AddEntry(h_all, "reconstructed pions", "f");
			legend->AddEntry(h_kaon, "mismatched kaons", "f");
			legend->AddEntry(h_other, "mismatched (other)", "f");
		}
		else{
			legend->AddEntry(h_all, "reconstructed kaons", "f");
			legend->AddEntry(h_kaon, "mismatched pions", "f");
			legend->AddEntry(h_other, "mismatched (other)", "f");
		}

		h_all->GetYaxis()->SetRangeUser(1,5e4);
		h_all->GetYaxis()->SetTitleOffset(1.11);

		c->SetLogy();


		h_all->Draw();
		h_other->Draw("SAME");
		h_kaon->Draw("SAME");

		legend->Draw();

		c->Print("Pion_kaon_mismatch_for_"+particle+".root");
		c->Print("Pion_kaon_mismatch_for_"+particle+".pdf");
		c->Print("Pion_kaon_mismatch_for_"+particle+".png");



}
